# -*- coding: utf-8 -*-
r"""
Product of formal power series

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
from itertools import imap
from sage.categories.formal_power_series import ExponentialPowerSeries, OrdinaryPowerSeries
from sage.combinat.species2.formal_power_series import FPS
from sage.combinat.species2.formal_power_series.operations.add import Add
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.arith import binomial
from sage.rings.integer import Integer


def clcall_private(ClassProd, cls, *args):

        # commutativity
        dic_fs = defaultdict(Integer)

        for i, (f, nf) in enumerate(args):
            # neutral element
            if f == f.category().one():
                continue
            # absorbing element
            elif f == f.category().zero():
                if nf != 0:
                    return f.category().zero()
            # distributivity f ⋅ (g + h) ⋅ r |--> f ⋅ g ⋅ r + f ⋅ h ⋅ r
            elif isinstance(f, Add):
                return Add(*tuple((ClassProd(*(tuple(dic_fs.items()) + ((g, ng),) + args[i+1:])), 1)
                                  for (g, ng) in f._dic_fs_.iteritems()), category=f.category())
            # associativity f ⋅ (g ⋅ h) |--> f ⋅ g ⋅ h
            elif isinstance(f, ClassProd):
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
                raise NotImplementedError("Inversion of fps I suppose... please implement...")

        if len(dic_fs.keys()) == 0:
            return f.category().one()
        elif len(dic_fs.keys()) == 1 and dic_fs.values()[0] == 1:
            return dic_fs.keys()[0]
        else:
            return super(ClassProd, cls).__classcall__(cls, tuple(dic_fs.items()))


class GenericProd(FPS):

    def is_pointing_of(self, g):
        """
        Test if `f` (*self*) is the pointing of `g`.

        In other terms, this method tests if

        MATH::

            f(t) = x \cdot g'(t)\,.

        """
        # TODO: this method is not consistent... should be implement in other operator class...
        from sage.combinat.species2.formal_power_series.operations.derivative import Derivative
        x = self.category().singletons()
        if len(self._dic_fs_) == 2 and x in self._dic_fs_.keys() and self._dic_fs_[x] == 1:
            dg = Derivative(g)
            return dg in self._dic_fs_.keys() and self._dic_fs_[dg] == 1
        return False

    def _valuation_(self):
        """
        Valuation of product of formal power series:

        MATH::

            val(f \times g) = val(f) + val(g)

        """
        return sum(map(lambda (f, nf): nf * f._valuation_(), self._dic_fs_.iteritems()))

    def _repr_(self):
        return "⋅".join(imap(lambda (ZF, nf): repr(ZF) + ("^%d" % nf if nf != 1 else ""), self._dic_fs_.iteritems()))


class ExponentialProd(GenericProd):
    """
    Product of exponential formal series

    MATH::

        (f \cdot g)(t) = f(t) \cdot g(t)

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args):
        return clcall_private(ExponentialProd, cls, *args)

    def __init__(self, dic_fs):
        FPS.__init__(self, category=ExponentialPowerSeries())
        self._dic_fs_ = dict(dic_fs)

    @cached_method
    def _fast_prod_pow_(self, f, k, n):
        """
        "Efficient" Computing of `[t^n]f(t)^k`.

        TESTS::

            sage: from sage.categories.formal_power_series import ExponentialPowerSeries
            sage: EPS = ExponentialPowerSeries()
            sage: x = EPS.singletons()
            sage: xx = x^5
            sage: xx._fast_prod_pow_(x, 5, 4)
            0
            sage: xx._fast_prod_pow_(x, 5, 5)
            120

        """
        if k == 0:
            res = Integer(1) if n == 0 else Integer(0)
        elif k % 2 == 0:
            res = sum(binomial(n, j) * self._fast_prod_pow_(f, k/2, j) * self._fast_prod_pow_(f, k/2, n-j)
                    for j in range(n+1))
        else:
            res = sum(binomial(n, j) * self._fast_prod_pow_(f, k-1, j) * f.coefficient(n-j) for j in range(n+1))
        return res

    @cached_method
    def coefficient(self, n):
        """
        (Naive) coefficient extraction

        TESTS::

            sage: from sage.categories.formal_power_series import ExponentialPowerSeries
            sage: EPS = ExponentialPowerSeries()
            sage: x = EPS.singletons()
            sage: o = EPS.one()
            sage: b = EPS.recursive_formal_power_series()
            sage: b.define(o + b*x*b)
            sage: list(b.coefficients(10)) # indirect doctest
            [1, 1, 4, 30, 336, 5040, 95040, 2162160, 57657600, 1764322560, 60949324800]
            sage: [catalan_number(n) * factorial(n) for n in range(11)]
            [1, 1, 4, 30, 336, 5040, 95040, 2162160, 57657600, 1764322560, 60949324800]

        """
        def rec_prod(fs, n):
            f = fs[0]
            if len(fs) == 1:
                return f.coefficient(n)

            acc = Integer(0)
            for i in range(f._valuation_(), n+1):
                fi = f.coefficient(i)
                if fi != 0:
                    acc += binomial(n, i) * fi * rec_prod(fs[1:], n-i)
            return acc

        fs = reduce(lambda e, f: e+f, tuple((f,)*nf for f, nf in self._dic_fs_.iteritems()), ())

        return rec_prod(fs, n)


class OrdinaryProd(GenericProd):

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args):
        return clcall_private(OrdinaryProd, cls, *args)

    def __init__(self, dic_fs):
        FPS.__init__(self, category=OrdinaryPowerSeries())
        self._dic_fs_ = dict(dic_fs)

    @cached_method
    def coefficient(self, n):
        """
        (Naive) coefficient extraction

        TESTS::

            sage: from sage.categories.formal_power_series import OrdinaryPowerSeries
            sage: OPS = OrdinaryPowerSeries()
            sage: x = OPS.singletons()
            sage: o = OPS.one()
            sage: b = OPS.recursive_formal_power_series()
            sage: b.define(o + b*x*b)
            sage: list(b.coefficients(10)) # indirect doctest
            [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]

        """
        def rec_prod(fs, n):
            f = fs[0]
            if len(fs) == 1:
                return f.coefficient(n)

            acc = Integer(0)
            for i in range(n+1):
                fi = f.coefficient(i)
                if fi != 0:
                    acc += fi * rec_prod(fs[1:], n-i)
            return acc

        fs = reduce(lambda e, f: e+f, tuple((f,)*nf for f, nf in self._dic_fs_.iteritems()), ())

        return rec_prod(fs, n)