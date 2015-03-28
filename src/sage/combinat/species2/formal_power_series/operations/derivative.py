# -*- coding: utf-8 -*-
r"""
Derivative of formal power series

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
from sage.combinat.species2.formal_power_series.operations.product import ExponentialProd, OrdinaryProd
from sage.combinat.species2.formal_power_series.operations.substitution import Substitution
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.integer import Integer


class Derivative(FPS):
    """
    Derivative of formal power series

    Properties:

     - `(f + g)'(t) = f'(t) + g'(t)`,
     - `(f \cdot g)'(t) = (f'\cdot g)(t) + (f \cdot g')(t)`,
     - `(f \circ g)'(t) = f'(g(t)) \cdot g'(t)`.

    EXAMPLES::

        sage: from sage.categories.formal_power_series import ExponentialPowerSeries
        sage: EPS = ExponentialPowerSeries()
        sage: o = EPS.one()
        sage: x = EPS.singletons()
        sage: b = EPS.recursive_formal_power_series()
        sage: b.define(o + b*x*b)
        sage: db = b.derivative()
        sage: list(db.coefficients(7))
        [1, 4, 30, 336, 5040, 95040, 2162160, 57657600]

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, f):
        # Z_E
        if f == f.category().sets():
            return f
        # product rule: (f ⋅ g)' --> f ⋅ g' + f' ⋅ g
        if isinstance(f, (ExponentialProd, OrdinaryProd)):
            Prod = f.__class__
            return Add(
                *[(Prod(*([(h, nh) if h != g else (h, nh-1)
                        for (h, nh) in f._dic_cis_.iteritems()] +
                          [(Derivative(g), Integer(1))])), ng)
                  for (g, ng) in f._dic_cis_.iteritems()]
            )
        # Additivity: (f + g)' --> f' + g'
        if isinstance(f, Add):
            return Add(*[(Derivative(ZG), ng) for (ZG, ng) in f._dic_cis_.iteritems()])
        # chain rule: (f(g))' --> f'(g) ⋅ g'
        if isinstance(f, Substitution):
            # f = g(h)
            g, h = f._f_, f._g_
            return Substitution(Derivative(g), h) * Derivative(h)
        # otherwise
        return super(Derivative, cls).__classcall__(cls, f)

    def __init__(self, f):
        FPS.__init__(self, category=f.category())
        self._f_ = f

    def _valuation_(self):
        return max(Integer(0), self._f_._valuation_() - 1)

    @lazy_attribute
    def coefficient(self):
        if self.category() == ExponentialPowerSeries():
            return self._coefficient_egs
        else:
            return self._coefficient_ogs

    def _coefficient_ogs(self, n):
        """
        TESTS::

            sage: from sage.categories.formal_power_series import ExponentialPowerSeries
            sage: OPS = OrdinaryPowerSeries()
            sage: o = OPS.one()
            sage: x = OPS.singletons()
            sage: b = OPS.recursive_formal_power_series()
            sage: b.define(o + b*x*b)
            sage: db = b.derivative()
            sage: list(db.coefficients(7))
            [1, 4, 15, 56, 210, 792, 3003, 11440]
            sage: [(n+1)*catalan_number(n+1) for n in range(8)]
            [1, 4, 15, 56, 210, 792, 3003, 11440]

        """
        return (n+1) * self._f_.coefficient(n+1)

    def _coefficient_egs(self, n):
        """
        TESTS::

            sage: from sage.categories.formal_power_series import ExponentialPowerSeries
            sage: EPS = ExponentialPowerSeries()
            sage: o = EPS.one()
            sage: x = EPS.singletons()
            sage: b = EPS.recursive_formal_power_series()
            sage: b.define(o + b*x*b)
            sage: db = b.derivative()
            sage: list(db.coefficients(7))
            [1, 4, 30, 336, 5040, 95040, 2162160, 57657600]
            sage: [factorial(n+1)*catalan_number(n+1) for n in range(8)]
            [1, 4, 30, 336, 5040, 95040, 2162160, 57657600]

        """
        return self._f_.coefficient(n+1)