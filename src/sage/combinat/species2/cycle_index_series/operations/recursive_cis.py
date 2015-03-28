# -*- coding: utf-8 -*-
"""
Cycle index series defined recursively

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
from sage.categories.formal_power_series import ExponentialPowerSeries, OrdinaryPowerSeries
from sage.misc.cachefunc import cached_method
from sage.rings.rational_field import QQ
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.combinat.sf.sf import SymmetricFunctions
from sage.structure.category_object import CategoryObject
from sage.structure.dynamic_class import dynamic_class


class RecursiveCIS(CategoryObject):
    """
    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: CIS = CycleIndexSeries()
        sage: ZZ = CIS.one()
        sage: ZX = CIS.singletons()
        sage: from sage.combinat.species2.cycle_index_series.operations.recursive_cis import RecursiveCIS
        sage: ZF = RecursiveCIS()
        sage: ZF.define(ZZ + ZF*ZX*ZF); ZF
        ZF := Z_1 + Z_X⋅ZF^2
        sage: [ZF.type_generating_series().coefficient(n) for n in range(10)]
        [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862]

        sage: TestSuite(ZF).run()

    """

    def __init__(self, name="ZF"):
        CategoryObject.__init__(self, category=CycleIndexSeries())
        self._def_ = None
        self._active_ = False
        self._name_ = name
        self._egs_ = ExponentialPowerSeries().recursive_formal_power_series(name=self._name_[1:])
        self._ogs_ = OrdinaryPowerSeries().recursive_formal_power_series(name=self._name_[1:])

        # ##### CUSTOM ####### used to inherit the parent methods
        base = self.__class__.__base__
        self.__class__ = dynamic_class("%s_with_category" % base.__name__,
                                       (self.__class__, self.category().parent_class, ),
                                       doccls=base)
        #####################

    def define(self, definition):
        """
        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: CIS = CycleIndexSeries()
            sage: ZZ = CIS.one()
            sage: ZX = CIS.singletons()
            sage: from sage.combinat.species2.cycle_index_series.operations.recursive_cis import RecursiveCIS
            sage: ZF = RecursiveCIS()
            sage: ZF.define(ZZ + ZF*ZX*ZF); ZF
            ZF := Z_1 + Z_X⋅ZF^2
            sage: ZF.define(ZF*ZX)
            Traceback (most recent call last):
            ...
            ValueError: This species is already defined: ZF := Z_1 + Z_X⋅ZF^2

        """
        if self._def_:
            raise ValueError("This cycle index series is already defined: %s" % repr(self))
        self._def_ = definition
        self._egs_.define(self._def_.generating_series())
        self._ogs_.define(self._def_.type_generating_series())
        self.Frobenius_characteristic = self._Frobenius_characteristic_def_

    def _repr_(self):
        if self._active_:
            return self._name_

        self._active_ = True
        s = self._name_ + " := " + repr(self._def_)
        self._active_ = False

        return s
    
    def Frobenius_characteristic(self, n):
        return SymmetricFunctions(QQ).zero()

    def _Frobenius_characteristic_def_(self, n):
        """

        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: CIS = CycleIndexSeries()
            sage: ZZ = CIS.one()
            sage: ZX = CIS.singletons()
            sage: from sage.combinat.species2.cycle_index_series.operations.recursive_cis import RecursiveCIS
            sage: ZB = RecursiveCIS(name="ZB")
            sage: ZB.define(ZZ + ZB*ZX*ZB); ZB # cycle index series of the binary tree species
            ZB := Z1 + ZX⋅ZB^2
            sage: s = SymmetricFunctions(QQ).s()
            sage: @cached_function
            ....: def bt(n):
            ....:   if n == 0: return s[[]]
            ....:   elif n == 1: return s[1]
            ....:   else: return sum(s[1] * bt(j) * bt(n-j-1) for j in range(n))
            sage: for n in range(6): assert(bt(n) == ZB.Frobenius_characteristic(n))

            sage: ZF = RecursiveCIS()
            sage: ZE = CIS.sets()
            sage: ZF.define(ZX*ZE.partitional_composite(ZF)) # cis of the rooted trees species
            sage: m = SymmetricFunctions(QQ).m()
            sage: for n in range(6): m(ZF.Frobenius_characteristic(n))
            0
            m[1]
            2*m[1, 1] + m[2]
            9*m[1, 1, 1] + 5*m[2, 1] + 2*m[3]
            64*m[1, 1, 1, 1] + 34*m[2, 1, 1] + 18*m[2, 2] + 13*m[3, 1] + 4*m[4]
            625*m[1, 1, 1, 1, 1] + 326*m[2, 1, 1, 1] + 171*m[2, 2, 1] + 119*m[3, 1, 1] + 63*m[3, 2] + 35*m[4, 1] + 9*m[5]

        """
        if self._valuation_() > n:
             return SymmetricFunctions(QQ).zero()

        if self._active_:
            if self._active_[-1] > n:
                self._active_.append(n)
                res = self._def_.Frobenius_characteristic(n)
                self._active_.pop()
            else:
                res = SymmetricFunctions(QQ).zero()
            return res

        self._active_ = [n]
        res = self._def_.Frobenius_characteristic(n)
        self._active_ = False

        return res

    def generating_series(self):
        return self._egs_

    @cached_method
    def type_generating_series(self):
        if not self._def_:
            raise AttributeError("This is not defined (yet)!")
        return self._ogs_