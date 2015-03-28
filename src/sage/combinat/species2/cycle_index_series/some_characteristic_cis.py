"""
Some characteristic cycle index series
"""
# *******************************************************************************
#       Copyright (C) 2015 Jean-Baptiste Priez <jbp@kerios.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
# *******************************************************************************
from sage.categories.formal_power_series import ExponentialPowerSeries, OrdinaryPowerSeries
from sage.combinat.partition import Partition
from sage.combinat.species2.cycle_index_series import CIS


class ZeroCIS(CIS):
    """
    The cycle index series `Z_0` of zero species:

    MATH::

        Z_0(p_1, p_2, p_3, \cdots) = 0

    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: Z0 = CycleIndexSeries().zero()
        sage: TestSuite(Z0).run()

    """

    def Frobenius_characteristic(self, n):
        """
        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: Z0 = CycleIndexSeries().zero()
            sage: sum(Z0.Frobenius_characteristic(n) for n in range(10))
            0

        """
        return self._sym_.zero()

    def _repr_(self):
        """
        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: CycleIndexSeries().zero()
            Z0

        """
        return "Z0"

    def generating_series(self):
        return ExponentialPowerSeries().zero()

    def type_generating_series(self):
        return OrdinaryPowerSeries().zero()


class OneCIS(CIS):
    """
    The cycle index series `Z_1` of one species:

    MATH::

        Z_1(p_1, p_2, p_3, \cdots) = 1

    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: Z1 = CycleIndexSeries().one()
        sage: TestSuite(Z1).run()

    """

    def Frobenius_characteristic(self, n):
        """
        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: Z1 = CycleIndexSeries().one()
            sage: sum(Z1.Frobenius_characteristic(n) for n in range(10))
            s[]

        """
        if n == 0:
            return self._sym_.one()
        return self._sym_.zero()

    def _repr_(self):
        """
        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: CycleIndexSeries().one()
            Z1

        """
        return "Z1"

    def generating_series(self):
        return ExponentialPowerSeries().one()

    def type_generating_series(self):
        return OrdinaryPowerSeries().one()


class SingletonsCIS(CIS):
    """
    The cycle index series `Z_X` of singletons species:

    MATH::

        Z_X(p_1, p_2, p_3, \cdots) = p_1

    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: ZX = CycleIndexSeries().singletons()
        sage: TestSuite(ZX).run()

    """

    def Frobenius_characteristic(self, n):
        """
        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: ZX = CycleIndexSeries().singletons()
            sage: sum(ZX.Frobenius_characteristic(n) for n in range(10))
            s[1]

        """
        if n == 1:
            return self._sym_.p().monomial(Partition([1]))
        return self._sym_.zero()

    def _repr_(self):
        """
        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: CycleIndexSeries().singletons()
            ZX

        """
        return "ZX"

    def generating_series(self):
        return ExponentialPowerSeries().singletons()

    def type_generating_series(self):
        return OrdinaryPowerSeries().singletons()


class SetsCIS(CIS):
    """
    The cycle index series `Z_E` of sets species:

    MATH::

        Z_E(p_1, p_2, p_3, \cdots) = \sum_{n \geqslant 0} h_n

    with `(h_n)` the complete symmetric functions basis.

    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: ZE = CycleIndexSeries().sets()
        sage: TestSuite(ZE).run()

    """

    def Frobenius_characteristic(self, n):
        """
        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: ZE = CycleIndexSeries().sets()
            sage: sum(ZE.Frobenius_characteristic(n) for n in range(10))
            h[] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + h[7] + h[8] + h[9]

        """
        return self._sym_.h().monomial(Partition([n]))

    def _repr_(self):
        """
        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: CycleIndexSeries().sets()
            ZE

        """
        return "ZE"

    def generating_series(self):
        return ExponentialPowerSeries().sets()

    def type_generating_series(self):
        return OrdinaryPowerSeries().sets()