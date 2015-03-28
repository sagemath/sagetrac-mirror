# -*- coding: utf-8 -*-
"""
Formal power series defined recursively

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
from sage.categories.formal_power_series import FormalPowerSeries
from sage.misc.cachefunc import cached_method
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.structure.dynamic_class import dynamic_class
from sage.structure.category_object import CategoryObject


class RecursiveFormalPowerSeries(CategoryObject):
    """
    Formal Power series defined recursively

    EXAMPLE::

        sage: from sage.categories.formal_power_series import ExponentialPowerSeries
        sage: EPS = ExponentialPowerSeries()
        sage: e = EPS.sets()
        sage: x = EPS.singletons()
        sage: r = EPS.recursive_formal_power_series(name="r")
        sage: r.define(x*e.substitution(r)) # rooted trees
        sage: list(r.coefficients(10))
        [0, 1, 2, 9, 64, 625, 7776, 117649, 2097152, 43046721, 1000000000]

        sage: f = EPS.recursive_formal_power_series(name="f")
        sage: f.define(e.substitution(x*f)) # forests
        sage: list(f.coefficients(10))
        [1, 1, 3, 16, 125, 1296, 16807, 262144, 4782969, 100000000, 2357947691]

        sage: from sage.categories.formal_power_series import OrdinaryPowerSeries
        sage: OPS = OrdinaryPowerSeries()
        sage: b = OPS.recursive_formal_power_series(name="b")
        sage: x = OPS.singletons()
        sage: o = OPS.one()
        sage: b.define(o + b*x*b) # binary trees
        sage: list(b.coefficients(10))
        [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]

    """

    def __init__(self, name="f", category=FormalPowerSeries()):
        CategoryObject.__init__(self, category=category)

        # ##### CUSTOM ####### used to inherit the parent methods
        base = self.__class__.__base__
        self.__class__ = dynamic_class("%s_with_category" % base.__name__,
                                       (self.__class__, self.category().parent_class, ),
                                       doccls=base)
        #####################

        self._def_ = None
        self._active_ = False
        self._name_ = name
        self._valuation = Infinity

    def define(self, definition):
        """

        """
        if self._def_:
            raise ValueError("This formal power series is already defined: %s" % repr(self))
        self._def_ = definition
        self._compute_valuation()

    def _compute_valuation(self):
        _old_valuation = None
        while self._valuation != _old_valuation:
            _old_valuation = self._valuation
            self._valuation = self._def_._valuation_()

    @cached_method
    def coefficient(self, n):
        """

        """
        if self._active_:
            if self._active_[-1] > n:
                self._active_.append(n)
                res = self._def_.coefficient(n)
                self._active_.pop()
            else:
                res = Integer(0)
            return res

        self._active_ = [n]
        res = self._def_.coefficient(n)
        self._active_ = False

        return res

    def _repr_(self):
        if self._active_:
            return self._name_

        self._active_ = True
        s = self._name_ + " := " + repr(self._def_)
        self._active_ = False

        return s

    def _valuation_(self):
        return self._valuation