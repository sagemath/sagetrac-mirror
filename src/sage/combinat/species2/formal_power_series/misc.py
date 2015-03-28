"""
Generic formal series associated to cycle index series
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
from sage.categories.formal_power_series import OrdinaryPowerSeries, ExponentialPowerSeries
from sage.combinat.partition import Partition
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.species2.formal_power_series import FPS
from sage.misc.cachefunc import cached_method
from sage.rings.arith import factorial, multinomial
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ


class genericOGS(FPS):
    """
    The specialization `\tilde{F}(x) = Z_F(x, x^2, x^3, \cdots)`.
    """

    def __init__(self, cis):
        FPS.__init__(self, category=OrdinaryPowerSeries())
        self._cis_ = cis

    def coefficient(self, n):
        z = self._cis_.Frobenius_characteristic(0)
        Q = z.base_ring()
        p = SymmetricFunctions(Q).p()
        h = SymmetricFunctions(Q).h()

        ch = self._cis_.Frobenius_characteristic(n)
        if ch.parent() == p:
            # ch(F[n]) is a sum of power sum
            res = ch.map_item(lambda I, c: (I.parent()([n]), c))
            if res == p.zero():
                return Integer(0)
            return res.coefficients()[0]
        # otherwise we convert into a homogeneous basis and h_lambda |--> 1
        res = h(ch).map_item(lambda I, c: (I.parent()([n]), c))
        if res == p.zero():
            return Integer(0)
        return res.coefficients()[0]

    @cached_method
    def _valuation_(self):
        n = 0
        zero = SymmetricFunctions(QQ).zero()
        # potentially infinite...
        while self._cis_.Frobenius_characteristic(n) == zero:
            n += 1
        return n


class genericEGS(FPS):
    """
    The specialization `F(x) = Z_F(x, 0, 0, \cdots)`.
    """

    def __init__(self, cis):
        FPS.__init__(self, category=ExponentialPowerSeries())
        self._cis_ = cis

    def coefficient(self, n):
        z = self._cis_.Frobenius_characteristic(0)
        Q = z.base_ring()
        p = SymmetricFunctions(Q).p()
        h = SymmetricFunctions(Q).h()

        def is_1k(I):
            return all(map(lambda i: i == 1, I))

        ch = self._cis_.Frobenius_characteristic(n)
        if ch.parent() == p:
            # ch(F[n]) is a sum of power sum
            res = ch.map_item(lambda I, c: (I, c*factorial(sum(I[:]))) if is_1k(I) else (I, 0))
            if res == p.zero():
                return Integer(0)
            return res.coefficients()[0]
        # otherwise we convert into a homogeneous basis and h_lambda |--> multinomial(lambda)
        res = h(ch).map_item(lambda I, c: (Partition([n]), c * multinomial(I[:])))
        if res == p.zero():
            return Integer(0)
        return res.coefficients()[0]

    @cached_method
    def _valuation_(self):
        n = 0
        zero = SymmetricFunctions(QQ).zero()
        # potentially infinite...
        while self._cis_.Frobenius_characteristic(n) == zero:
            n += 1
        return n