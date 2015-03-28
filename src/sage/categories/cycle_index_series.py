# -*- coding: utf-8 -*-
r"""
Cycle index series

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
from sage.categories.objects import Objects
from sage.rings.infinity import Infinity
from sage.categories.category import Category
from sage.combinat.sf.sf import SymmetricFunctions
from sage.misc.abstract_method import abstract_method


class CycleIndexSeries(Category):
    """
    The *cycle index series* of a species of structures `F` is the formal power
    series (in the symmetric functions)

    MATH::

        Z_F(p_1, p_2, p_3, \cdots) = \sum_{n \geqslant 0} \frac{1}{n!} \left(
            \sum_{\sigma \in \mathcal{S}_n} \mathtt{fix}\, F[\sigma]\, p_1^{\sigma_1} p_2^{\sigma_2} p_3^{\sigma_3}
            \cdots
        \right)

    where `\mathcal{S}_n` denotes the group of permutations of `[n]` and `\mathtt{fix}\, F[\sigma]` is the number of
    `F`-structures on `[n]` fixed by `F[\sigma]`.

    (section 1.2, Definition 6, _[BBL])
    """

    def super_categories(self):
        return [Objects()]

    def zero(self):
        from sage.combinat.species2.cycle_index_series.some_characteristic_cis import ZeroCIS
        return ZeroCIS()

    def one(self):
        from sage.combinat.species2.cycle_index_series.some_characteristic_cis import OneCIS
        return OneCIS()

    def singletons(self):
        from sage.combinat.species2.cycle_index_series.some_characteristic_cis import SingletonsCIS
        return SingletonsCIS()

    def sets(self):
        from sage.combinat.species2.cycle_index_series.some_characteristic_cis import SetsCIS
        return SetsCIS()

    class ParentMethods:

        def _test_cycle_index_series(self, **options):
            tester = self._tester(**options)
            z = self.Frobenius_characteristic(0)
            Q = z.base_ring()
            tester.assertIn(z, SymmetricFunctions(Q))

        @abstract_method
        def Frobenius_characteristic(self, n):
            """
            The *Frobenius characteristic* of the natural symmetric group action on `F[n]`.

            MATH::

                ch(F[n]) = [n]Z_F(p_1, p_2, p_3, \cdots) =
                \sum_{\sigma \in \mathcal{S}_n} \mathtt{fix}\, F[\sigma]\, p_1^{\sigma_1} p_2^{\sigma_2} p_3^{\sigma_3}

            :param n: an non-negative integer
            """

        def generating_series(self):
            """
            The generating series of the species of structures `F` is the formal power series

            MATH::

                F(t) = \sum_{n \geqslant 0} f_n \frac{t^n}{n!}\,,

            where `f_n = |F[n]|` the number of `F`-structures on a set of cardinality `n`.

            This generating series is a specialization of the *cycle index series*:

            MATH::

                F(t) = Z_F(t, 0, 0, \cdots)\,.

            (section 1.2, Theorem 8, _[BBL])
            """
            from sage.combinat.species2.formal_power_series.misc import genericEGS
            return genericEGS(self)

        def exponential_generating_series(self): return self.generating_series()
        def egs(self): self.generating_series()

        def type_generating_series(self):
            """
            The (isomorphism) type generating series of the species of structures `F` is the formal power series

            MATH::

                \tilde{F}(t) = \sum_{n \geqslant 0} \tilde{f_n} t^n\,,

            where `\tilde{f_n} = |T(F[n])|` the number of equivalence classes associated to `F`-structures on a set of
            cardinality `n`.

            This generating series is a specialization of the *cycle index series*:

            MATH::

                F(t) = Z_F(t, t^2, t^3, \cdots)\,.

            (section 1.2, Theorem 8, _[BBL])
            """
            from sage.combinat.species2.formal_power_series.misc import genericOGS
            return genericOGS(self)

        def isomorphism_type_generating_series(self): return self.type_generating_series()
        def ordinary_generating_series(self): return self.type_generating_series()
        def ogs(self): return self.type_generating_series()

        def add(ZF, ZG):
            """
            The *sum of cycle index series* `Z_F` and `Z_G`:

            MATH::

                (Z_F + Z_G)(p_1, p_2, p_3, \cdots) = Z_F(p_1, p_2, p_3, \cdots) + Z_G(p_1, p_2, p_3, \cdots)

            (section 1.3, _[BBL])

            :param ZF, ZG: are cycle index series
            """
            from sage.combinat.species2.cycle_index_series.operations.add import Add
            return Add((ZF,1), (ZG,1))

        __add__ = add

        def product(ZF, ZG):
            """
            The *product of cycles index series* `Z_F` and `Z_G`:

            MATH::

                (Z_F \cdot Z_G)(p_1, p_2, p_3, \cdots) = Z_F(p_1, p_2, p_3, \cdots) \cdots Z_G(p_1, p_2, p_3, \cdots)

            (section 1.3, _[BBL])

            :param ZF, ZG: are cycle index series
            """
            from sage.combinat.species2.cycle_index_series.operations.product import Prod
            return Prod((ZF, 1), (ZG, 1))

        __mul__ = _mul_ = product


        def restriction(self, min=0, max=Infinity):
            """
            The *restriction* of the cycle index series `Z_F` to order `n` with `min \leqslant n \leqslant max`.

            :param min and max: non-negative integers
            """
            from sage.combinat.species2.cycle_index_series.operations.restriction import RestrictionCIS
            return RestrictionCIS(self, min, max)

        def partitional_composite(ZF, ZG):
            """
            The *(partitional) composite* of `Z_G` into `Z_F`.

            :param ZF and ZG: cycles index series

            """
            from sage.combinat.species2.cycle_index_series.operations.partitional_composite import Composite
            return Composite(ZF, ZG)

        __call__ = partitional_composite

        composite = partitional_composite

        def derivative(ZF):
            """
            The *derivative* of `Z_F`.

            MATH::

                Z_F'(p_1, p_2, p_3, \cdots) = \left(\frac{\partial}{\partial p_1} Z_F\right)(p_1, p_2, p_3, \cdots)\,.

            """
            from sage.combinat.species2.cycle_index_series.operations.derivative import Derivative
            return Derivative(ZF)

        def pointing(ZF):
            """
            The *pointing* of `Z_F`:

            MATH::

                Z_F^\bullet(p_1, p_2, p_3, \cdots) = p_1 \left(
                    \frac{\partial}{\partial p_1} Z_F
                \right) (p_1, p_2, p_3, \cdots)\,.

            """
            return ZF.category().singletons() * ZF.derivative()

        def Hadamard_product(ZF, ZG):
            """
            The *Hadamard product* of cycles index series

            MATH::

                (Z_F \times Z_G)(p_1, p_2, p_3, \cdots) = Z_F(p_1, p_2, p_3, \cdots) \times Z_G(p_1, p_2, p_3, \cdots)\,.

            """
            from sage.combinat.species2.cycle_index_series.operations.hadamard_product import HadamardProduct
            return HadamardProduct((ZF, 1), (ZG, 1))

        hadamard_product = Hadamard_product


        def functorial_composite(self, ZG):
            """
            The *functorial composite* of cycle index series

            MATH::

                (Z_F \Box Z_G)(p_1, p_2, p_3, \cdots) = ??

            """
            from sage.combinat.species2.cycle_index_series.operations.functorial_composite import FunctorialComposite
            return FunctorialComposite(self, ZG)

        __getitem__ = functorial_composite

        def arithmetic_product(ZF, ZG):
            """

            Reference:
            ----------
            .. On the arithmetic product of combinatorial species
               Manuel Maia and Miguel Méndez
               2008

            """
            # TODO
            NotImplemented

        def _valuation_(self):
            return self.generating_series()._valuation_()
