# -*- coding: utf-8 -*-
r"""
Hadamard product of formal power series

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
from collections import defaultdict
from sage.categories.formal_power_series import FormalPowerSeries
from sage.combinat.species2.formal_power_series import FPS
from sage.combinat.species2.formal_power_series.operations.add import Add
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.integer import Integer


class HadamardProduct(FPS):
    """
    Hadamard product of formal power series

    MATH::

        (f \times g)(t) = \sum_{\geqslant 0} f_n g_n t^n\,.

    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):

        # commutativity
        dic_fs = defaultdict(Integer)

        for i, (f, nf) in enumerate(args):
            # neutral element
            if f == f.category().sets():
                continue
            # absorbing element
            elif f == f.category().zero():
                if nf != 0:
                    return f.category().zero()
                # else continue
            # distributive
            elif isinstance(f, Add):
                return Add(*tuple((HadamardProduct(*(tuple(dic_fs.items()) + ((g, ng),) + args[i+1:])), 1)
                                  for (g, ng) in f._dic_fs_.iteritems()))
            # associative
            elif isinstance(f, HadamardProduct):
                for (g, ng) in f._dic_fs_.iteritems():
                    dic_fs[g] += ng
            # otherwise
            else:
                dic_fs[f] += nf

        # # cleaning ##
        for (f, nf) in list(dic_fs.items()):
            if nf == 0:
                del dic_fs[f]
            if nf < 0:
                raise NotImplementedError("Virtual species I suppose... please implement...")

        if len(dic_fs.keys()) == 0:
            return args[0][0].category().sets()
        elif len(dic_fs.keys()) == 1 and dic_fs.values()[0] == 1:
            return dic_fs.keys()[0]
        else:
            return super(HadamardProduct, cls).__classcall__(cls, tuple(dic_fs.items()), **opts)

    def __init__(self, dic_fs, category=FormalPowerSeries()):
        FPS.__init__(self, category=category)
        self._dic_fs_ = dict(dic_fs)

    def _valuation_(self):
        return max(map(lambda (f, nf): f._valuation_(), self._dic_fs_.iteritems()))

    def coefficient(self, n):
        if n < self._valuation_():
            return 0
        return reduce(lambda acc, (f, nf): acc * f.coefficient(n) ** nf, self._dic_fs_.iteritems(), Integer(1))