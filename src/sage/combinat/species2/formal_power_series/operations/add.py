# -*- coding: utf-8 -*-
r"""
Sum of formal power series

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998
}
"""
# *****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************
from collections import defaultdict
from itertools import imap
from sage.categories.formal_power_series import FormalPowerSeries
from sage.combinat.species2.formal_power_series import FPS
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.integer import Integer


class Add(FPS):
    """
    Sum of power series

    MATH::

        (f + g)(t) = f(t) + g(t)\,.

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        cat = opts["category"]
        # args = ((f, nf), (g, ng), ...) means nf.f + ng.g + ...
        dic_fs = defaultdict(Integer)
        for (f, nf) in args:
            # Neutral element:
            if f == cat.zero():
                continue
            # Associativity
            elif isinstance(f, Add):
                for (g, ng) in f._dic_series_.iteritems():
                    dic_fs[g] += ng * nf
            else:
                dic_fs[f] += nf

        # ### simplify (difference of fps) ####
        for (f, nf) in list(dic_fs.items()):
            if nf == 0:
                del dic_fs[f]
        ######################################

        if len(dic_fs.keys()) == 0:
            return cat.zero()
        elif len(dic_fs.keys()) == 1 and dic_fs[dic_fs.keys()[0]] == 1:
            return dic_fs.keys()[0]
        else:
            return super(Add, cls).__classcall__(cls, tuple(dic_fs.items()), **opts)

    def __init__(self, dic_fs, category=FormalPowerSeries()):
        FPS.__init__(self, category=category)
        self._dic_fs_ = dict(dic_fs)

    def coefficient(self, n):
        """
        MATH::

            [t^n](f + g)(t) = [t^n]f(t) + [t^n]g(t)

        :param n: an integer
        :return: `[t^n](f + g)(t)`.
        """
        if n < self._valuation_():
            return Integer(0)
        return sum(imap(lambda (f, nf): nf * f.coefficient(n),
                        self._dic_fs_.iteritems()))

    def _repr_(self):
        return " + ".join(imap(lambda (f, nf): (repr(nf) + "⋅" if nf != 1 else "") + repr(f),
                               self._dic_fs_.iteritems()))

    def _valuation_(self):  # TODO: cached or not cached?
        return min(*map(lambda f: f._valuation_(), self._dic_fs_.keys()))