# -*- coding: utf-8 -*-
"""
Cycle index series defined recursively

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998
"""
#*****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.rings.rational_field import QQ
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.combinat.sf.sf import SymmetricFunctions
from sage.structure.parent import Parent


class RecursiveCIS(Parent):
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

    """

    def __init__(self):
        Parent.__init__(self, category=CycleIndexSeries())
        self._def_ = None
        self._active_ = False

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
            raise ValueError("This species is already defined: %s"%repr(self))
        self._def_ = definition

    def _repr_(self):
        if self._active_:
            return "ZF"

        self._active_ = True
        s = "ZF := " + repr(self._def_)
        self._active_ = False

        return s

    def Frobenius_characteristic(self, n):
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
            sage: [ZF.type_generating_series().coefficient(n) for n in range(10)] #indirect doctest
            [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862]
            sage: [ZF.generating_series().coefficient(n)/factorial(n) for n in range(10)] #indirect doctest
            [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862]

        """
        if self._active_:
            if self._active_[-1] > n:
                self._active_.append(n)
                res = self._def_.Frobenius_characteristic(n)
                self._active_.pop()
            else:
                res =  SymmetricFunctions(QQ).zero()
            return res

        self._active_ = [n]
        res = self._def_.Frobenius_characteristic(n)
        self._active_ = False

        return res

