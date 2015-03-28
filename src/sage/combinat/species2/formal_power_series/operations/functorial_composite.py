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
from sage.combinat.species2.formal_power_series import FPS
from sage.combinat.species2.formal_power_series.operations.add import Add
from sage.combinat.species2.formal_power_series.operations.hadamard_product import HadamardProduct
from sage.combinat.species2.formal_power_series.operations.product import ExponentialProd
from sage.misc.classcall_metaclass import ClasscallMetaclass


class FunctorialComposite(FPS):
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

    def _valuation_(self):
        return self._g_.coefficient(self._g_._valuation_())

    def coefficient(self, n):
        """
        MATH::

            [t^n] (f \Box g)(t) = f_{g_n}\,.

        TESTS::

            sage:
        """
        return self._f_.coefficient(self._g_.coefficient(n))