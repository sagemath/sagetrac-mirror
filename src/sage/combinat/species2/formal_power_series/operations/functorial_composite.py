# -*- coding: utf-8 -*-
r"""
Functorial composite of formal power series

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
from sage.combinat.species2.formal_power_series import FPS, ValuationFPS
from sage.combinat.species2.formal_power_series.operations.add import Add
from sage.combinat.species2.formal_power_series.operations.hadamard_product import HadamardProduct
from sage.combinat.species2.formal_power_series.operations.product import ExponentialProd
from sage.misc.classcall_metaclass import ClasscallMetaclass


class FunctorialComposite(ValuationFPS, FPS):
    """
    Functorial composite of formal power series

    MATH::

        (f \Box g)(t) = \sum_{n \seq 0} f_{g_n} t^n\,.

    Properties:

     - `(f \Box (g \Box h))(t) = ((f \Box g) \Box h)(t)`,
     - neutral element `e^\bullet (= x e')`,
     - Distributive on the right by the operation of hadamard product:
     `((f \times g) \Box h)(t) = ((f \Box h) \times (g \Box h))(t)`.

    EXAMPLES::


    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, f, g):

        # Neutral element: E• = X⋅E'
        if isinstance(f, ExponentialProd) and f.is_pointing_of(ExponentialPowerSeries().sets()):
            return g
        if isinstance(g, ExponentialProd) and g.is_pointing_of(ExponentialPowerSeries().sets()):
            return f
        # Associative
        if isinstance(g, FunctorialComposite):
            return FunctorialComposite(FunctorialComposite(f, g._f_), g._g_)
        # distributive on the right: x
        if isinstance(f, HadamardProduct):
            return HadamardProduct(*tuple((FunctorialComposite(h, g), nh)
                                          for (h, nh) in f._dic_fs_.iteritems()))
        # distributive on the right: +
        if isinstance(f, Add):
            return Add(*tuple((FunctorialComposite(h, g), nh)
                              for (h, nh) in f._dic_fs_.iteritems()))
        # otherwise
        return super(FunctorialComposite, cls).__classcall__(cls, f, g)

    def __init__(self, f, g):
        FPS.__init__(self, category=ExponentialPowerSeries())
        self._f_, self._g_ = f, g

        ValuationFPS.__init__(self)
        self._valuation_registration_([f, g])
        self._valuation_update_()

    def _valuation_compute_(self):
        # FIXME min {n | [t^n]g(t) >= val(f)} ??
        f, g = self._f_, self._g_
        n = 0
        while g.coefficient(n) < f.valuation():
            n += 1
        return n

    def coefficient(self, n):
        """
        MATH::

            [t^n] (f \Box g)(t) = f_{g_n}\,.

        TESTS::

            Sp = Species()
            sage: E  = Sp.sets()
            sage: S = E*E
            sage: E2 = E.restricted(min=2, max=2)
            sage: P2 = E2*E
            sage: G = S.functorial_composite(P2)
            sage: list(G.egs().coefficients(5))
            [1, 1, 2, 8, 64, 1024]

        """
        return self._f_.coefficient(self._g_.coefficient(n))